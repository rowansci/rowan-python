"""Tests for the user model."""

from rowan.user import User


def test_credits_available_string_shows_organization_pool() -> None:
    """Users have no weekly credits of their own; only the organization's are shown."""
    user = User.model_validate(
        {
            "uuid": "user-uuid",
            "username": "ada",
            "email": "ada@example.com",
            "credits": 0,
            "weekly_credits": 20,
            "organization": {"name": "Lab", "weekly_credits": 100, "credits": 5},
        }
    )

    assert not hasattr(user, "weekly_credits")
    assert user.credits_available_string() == (
        "Credits: 0.0\nOrganization Weekly Credits: 100.0\nOrganization Credits: 5.0"
    )
